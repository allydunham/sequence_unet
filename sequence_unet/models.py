"""
Load, download and initialise trained and untrained Sequence UNET models.
"""
from logging import error
import os
import ftplib
import urllib.request
import tqdm
from shutil import unpack_archive
from contextlib import closing
from tensorflow.keras import layers, models

from sequence_unet.graph_cnn import GraphCNN
from sequence_unet import metrics

__all__ = ["sequence_unet", "cnn_top_model", "download_trained_model",
           "download_all_models", "load_trained_model", "MODELS", "CUSTOM_OBJECTS",
           "BIOSTUDIES_FTP"]

CUSTOM_OBJECTS = {
    "masked_binary_crossentropy": metrics.masked_binary_crossentropy,
    "masked_accuracy": metrics.masked_accuracy,
    "WeightedMaskedBinaryCrossEntropy": metrics.WeightedMaskedBinaryCrossEntropy,
    "GraphCNN": GraphCNN
}
"""
Dictionary containing the object mappings required to load Sequence UNET
Keras models using `tf.keras.models.load_model`.
"""

BIOSTUDIES_FTP = "biostudies/nfs/S-BSST/732/S-BSST732"
"""
FTP path to the BioStudies directory containing the model data
"""

BIOSTUDIES_HTTPS = "https://www.ebi.ac.uk/biostudies/files/S-BSST732"
"""
Base HTTPS path for models on Biostudies
"""

MODELS = [
    'freq_classifier',
    'pregraph_freq_classifier',
    'pssm_predictor',
    'pregraph_pssm_predictor',
    'patho_top',
    'pregraph_patho_top',
    'patho_finetune',
    'pregraph_patho_finetune'
]
"""
List of model IDs for each trained Sequence UNET model. HDF5 (.h5) and SavedModel (.tf.tar.gz) files are available for each model at {BIOSTUDIES_FTP}/Files/{name}.{ext}.
"""

def _make_progess_file_writer(file, pbar):
    def f(data):
        file.write(data)
        pbar.update(len(data))
    return f

def _download_model_ftp(model, ext, dl_path):
    """
    Download a model from BioStudies over FTP
    """
    with closing(ftplib.FTP("ftp.ebi.ac.uk")) as ftp:
        ftp.login()
        ftp.cwd(f"{BIOSTUDIES_FTP}/Files/")
        file_size = ftp.size(f"{model}.{ext}")
        pbar = tqdm.tqdm(desc=f"Downloading {model}", total=file_size, unit="B", unit_scale=True)
        with open(dl_path, "w+b") as dl_file, pbar as pbar:
            writer = _make_progess_file_writer(dl_file, pbar)
            res = ftp.retrbinary(f"RETR {model}.{ext}", writer)

        if not res.startswith('226 Transfer complete.'):
            raise ConnectionError(f"FTP download failed with response \"{res}\"")

def _download_model_http(model, ext, dl_path):
    """
    Download a model from BioStudies over HTTP
    """
    # Based on https://stackoverflow.com/a/41107237
    url = f"{BIOSTUDIES_HTTPS}/{model}.{ext}"
    with closing(urllib.request.urlopen(url)) as res:
        if not res.status == 200:
            raise ConnectionError(f"HTTP connection failed with status code {res.status}")

        file_size = res.getheader('content-length')
        file_size = int(file_size) if file_size else 550000000
        blocksize = max(4096, file_size//100)
        pbar = tqdm.tqdm(desc=f"Downloading {model}", total=file_size, unit="B", unit_scale=True)
        with open(dl_path, "w+b") as dl_file, pbar as pbar:
            writer = _make_progess_file_writer(dl_file, pbar)
            while True:
                chunk = res.read(blocksize)
                if not chunk:
                    break
                writer(chunk)

def download_trained_model(model, root=".", model_format="tf", use_ftp=True):
    """
    Download a trained Sequence UNET model.

    Download the specified trained Sequence UNET model from BioStudies.
    Models are specified using the IDs indicated in models.MODELS, which also maps them to BioStudies files.
    Each model comes in a sequence only version (X) and version accepting additional structural input (pregraph_X).
    FTP download is recommended to give the best transfer speeds and lowest load for
    BioStudies, but an HTTP fallback is also provided if FTP isn't possible.

    Available models:

    :freq_classifier: Classifier predicting where variants occur above or below 0.01 observation frequency in a cross species multiple sequence alignment, as a proxy for deleteriousness.
    :pssm_predictor: Model predicting multiple alignment frequencies, which can be converted into a PSSM.
    :patho_top: Classifier predicting variant pathogenicity, trained as a new classifier head for ``freq_classifier`` on ClinVar data.
    :patho_finetune: Classifier predicting variant pathogenicity, trained by finetuning `freq_classifier` on ClinVar data.
    :pregraph_freq_classifier: Equivalent to ``freq_classifier`` taking structural input.
    :pregraph_pssm_predictor: Equivalent to ``pssm_predictor`` taking structural input.
    :pregraph_patho_top: Equivalent to ``patho_top`` taking structural input.
    :pregraph_patho_finetune: Equivalent to ``patho_finetune`` taking structural input.

    Parameters
    ----------
    model        : str
        Sequence UNET model to download (see options in description and `MODELS`).
    root         : str
        Root directory to download to.
    model_format : {'tf', 'h5'}
        Format to download the models in.
    use_ftp      : bool
        Download over FTP rather than HTTP.

    Returns
    -------
    str
        Path to downloaded model.
    """
    if not os.path.isdir(root):
        raise FileNotFoundError(f"No directory found at {root}")

    if model not in MODELS:
        raise ValueError(f"model not recognised. Must one of: {', '.join(MODELS)}")

    if not model_format in ('tf', 'h5'):
        raise ValueError(f"Model format ({model_format}) not recognised. Must be tf or h5")

    ext = "tf.tar.gz" if model_format == 'tf' else 'h5'
    dl_path = f"{root}/{model}.{ext}"

    dl_func = _download_model_ftp if use_ftp else _download_model_http
    try:
        dl_func(model, ext, dl_path)
    except (Exception, KeyboardInterrupt) as err:
        try:
            os.remove(dl_path)
        except FileNotFoundError:
            pass
        except Exception as err:
            error((f"Failed to delete partial download with error {err}. "
                    f"There may be an incomplete download left at {dl_path}."))
        raise ConnectionError((f"{'FTP' if use_ftp else 'HTTP'} download failed for {model}. "
                               f"Consider trying {'HTTP' if use_ftp else 'FTP'} with "
                               f"use_ftp={not use_ftp}"))

    if model_format == "tf":
        unpack_archive(dl_path, format="gztar")
        os.remove(dl_path) # Remove Tar archive
        dl_path = dl_path[:-7] # Drop .tar.gz extension

    return dl_path

def download_all_models(root=".", model_format="tf", use_ftp=True):
    """
    Download all trained Sequence UNET models.

    Download all trained Sequence UNET models from BioStudies. See `download_trained_model` for a description of the available models.
    FTP download is recommended to give the best transfer speeds and lowest load for
    BioStudies, but an HTTP fallback is also provided if FTP isn't possible.

    Parameters
    ----------
    root         : str
        Root directory to download to.
    model_format : {'tf', 'h5'}
        Format to download the models in.
    use_ftp      : bool
        Download over FTP rather than HTTP.

    Returns
    -------
    Dict
        Dictionary of paths to the downloaded models, keyed by the model
    """
    if not os.path.isdir(root):
        raise FileNotFoundError(f"No directory found at {root}")

    if not model_format in ('tf', 'h5'):
        raise ValueError(f"Model format ({model_format}) not recognised. Must be tf or h5")

    paths = {}
    for model in MODELS:
        paths[model] = download_trained_model(model, root=root, model_format=model_format,
                                              use_ftp=use_ftp)

    return paths

def load_trained_model(model, root=".", download=False, model_format="tf", use_ftp=True):
    """
    Load Sequence UNET models.

    Load a trained Sequence UNET model downloaded with `download_trained_model` or directly from BioStudies. This function provides a convenient wrapper around the Keras load_model function, allowing path or model name input and passing the appropriate custom_objects dictionary.
    If using load_model directly the `CUSTOM_OBJECTS` dictionary is required so TensorFlow can locate the custom Sequence UNET layers and metrics.

    Parameters
    ----------
    model        : str
        Model to load. This can be a direct path or a model name (see `MODELS`), with paths taking precedence. If a name is passed it will be searched for in root.
    root         : str
        Root directory to locate models passed by name. It is ignored if a full path is passed.
    download     : bool
        Download the requested model if it is not located.
    model_format : str
        Format for the model if model is an ID rather than a full path.
    use_ftp      : bool
        Download over FTP rather than HTTP.

    Returns
    -------
    tf.keras.Model
        The trained Sequence UNET model
    """
    if os.path.exists(model):
        path = model

    elif model in MODELS:
        path = f"{root}/{model}.{model_format}"
        if not os.path.exists(path):
            if download:
                download_trained_model(model, root, model_format=model_format, use_ftp=use_ftp)
            else:
                raise FileNotFoundError("No model found at {path}")

    else:
        raise ValueError(("Unrecognised model - must be a path to a SavedModel or one of the "
                          "options in sequence_unet.models.MODELS"))

    return models.load_model(path, custom_objects=CUSTOM_OBJECTS)

def sequence_unet(filters=8, kernel_size=5, num_layers=4, output_size=20,
                  dropout=0, graph_layers=None, graph_activation="relu",
                  conv_activation="relu", pred_activation="sigmoid",
                  kernel_regulariser=None, batch_normalisation=False):
    """
    Initialise a Sequence UNET TensorFlow Keras model.

    Generate a functional style Keras Sequence UNET model.
    This model takes one-hot encoded sequence and optionally inverted residue disance matrix structural information (see `graph_cnn.contact_graph`) and generates a matrix of predictions for each position in the sequence.
    It is a 1D CNN based model that uses a U-shaped compression and decompression architecture to spread information through the sequence.

    Parameters
    ----------
    filters              : int
        Number of convolutional filters on the top layers. Lower layers have F x 2^N filters, so large numbers of filters in deep networks quickly scales to very many weights.
    kernel_size          : int
        Width of 1D convolutional kernals.
    num_layers           : int
        Number of down/up sampling layers.
    output_size          : int
        Number of output predictions for each position in the sequence. I.e. the number of columns in the output matrix.
    dropout              : float
        Proportion of neurons dropped out between layers.
    graph_layers         : int or None
        Number of GraphCNN layers preprocessing structural features to feed into the main network.
    graph_activation     : str
        Activation function for GraphCNN layers.
    conv_activation      : str
        Activation function for CNN layers.
    pred_activation      : str
        Activation function for final layer. Sigmoid was used for frequency classification and softmax for PSSM frequency prediction.
    kernel_regulariser   : str or None
        Type of kernel regularisation to apply (l1, l2, l1/l2). See Keras Conv1D kernel_regularizer options for details.
    batch_normalisation  : bool
        Apply batch normalisation between contraction and upsampling layers.

    Returns
    -------
    tf.keras.Model
        The Sequence UNET model

    Notes
    -----
    Expected input:
        (B, N, 20) arrays with B batches of Nx20 matrices, representing
        sequences with one hot encoding on each row. An additional (B, N, N)
        array with NxN contact graphs is required when graph_layers is not none.
        The UNET architecture means N must be divisble by 2^(num_layers - 1).
        This can be achieved with padding and masking.

    Output:
        (B, N, 20) arrays of B batches of Nx20 arrays, with
        the predicted features of each AA at each position as rows.

    """
    num_layers = num_layers - 1 # 0 index layers

    input_seq = layers.Input(shape=[None, 20], name='input_seq')
    inputs = input_seq
    x = input_seq

    if graph_layers is not None:
        input_graph = layers.Input(shape=[None, None], name='input_graph')
        inputs = [input_seq, input_graph]
        for i, units in enumerate(graph_layers):
            x = GraphCNN(units, activation=graph_activation, name=f'graph_{i}')([input_graph, x])

        x = layers.concatenate([input_seq, x], name="graph_concat")

    # Contraction
    contraction = []
    for i in range(num_layers):
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'down_{i}_conv_1',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'down_{i}_conv_2',
                          kernel_regularizer=kernel_regulariser)(x)
        contraction.append(x)
        x = layers.MaxPool1D(pool_size=2, strides=2, padding='same', name=f'down_{i}_max_pool')(x)

        if batch_normalisation:
            x = layers.BatchNormalization(name=f"down_{i}_batch_normalisation")(x)

        if dropout:
            x = layers.SpatialDropout1D(dropout, name=f'down_{i}_dropout')(x)

    # Bottom layer
    x = layers.Conv1D(filters * 2 ** num_layers, kernel_size, 1,
                      padding='same', activation=conv_activation, name=f'bottom_conv_1',
                      kernel_regularizer=kernel_regulariser)(x)
    x = layers.Conv1D(filters * 2 ** num_layers, kernel_size, 1,
                      padding='same', activation=conv_activation, name=f'bottom_conv_2',
                      kernel_regularizer=kernel_regulariser)(x)

    if batch_normalisation:
            x = layers.BatchNormalization(name="bottom_batch_normalisation")(x)

    if dropout:
        x = layers.SpatialDropout1D(dropout, name='bottom_dropout')(x)

    # Expansion
    for i in range(num_layers - 1, -1, -1):
        x = layers.UpSampling1D(2, name=f'up_{i}_upsample')(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_1',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.concatenate([contraction[i], x], name=f"up_{i}_concat")

        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_2',
                          kernel_regularizer=kernel_regulariser)(x)
        x = layers.Conv1D(filters * 2 ** i, kernel_size, 1, padding='same',
                          activation=conv_activation, name=f'up_{i}_conv_3',
                          kernel_regularizer=kernel_regulariser)(x)

        if batch_normalisation:
            x = layers.BatchNormalization(name=f"up_{i}_batch_normalisation")(x)

        if dropout:
            x = layers.SpatialDropout1D(dropout, name=f'up_{i}_dropout')(x)

    preds = layers.Conv1D(output_size, 1, 1, activation=pred_activation, name="predictor",
                          kernel_regularizer=kernel_regulariser)(x)

    return models.Model(inputs=inputs, outputs=preds)

def cnn_top_model(bottom_model, features=True, output_size=20,
                  tune_layers=3, kernel_size=3,
                  activation="sigmoid", kernel_regulariser=None,
                  activity_regulariser=None, dropout=0):
    """
    Add a 1D CNN top model to a trained Tensorflow Keras model

    Add a fresh prediction head onto a trained model to finetune it for a new application. This
    was used to create the Sequence UNET ClinVar pathogenicity prediction top model. The new layer is a 1D CNN designed to interface with Sequence UNET architechture. A similar approach can easily be used to add a different final classification head.

    bottom_model: PSSM model to train on top of. Bottom_model=None creates a simple CNN model.
    features: Start from model features rather than predicted output.
    tune_layers: Number of top layers to set as trainable.

    Parameters
    ----------
    bottom_model          : tf.keras.Model
        Model to train on top of.
    features              : bool
        Use the final feature layer as input to the top model, replacing the current prediction layer, instead of starting from the current prediction layer.
    output_size           : int
        Number of filters in the new output layer. This is the number of predictions made for each sequence position.
    tune_layers           : int
        Number of layers in the bottom model to unfreeze and train. Use 0 to only train the new layer and -1 to keep all layers trainable.
    kernel_size           : int
        CNN kernel width
    activation            : str
        Activation function for the top layer, determining the type of prediction made.
    kernel_regulariser    : str or None
        Kernel regularisation to apply (e.g. l1, l2, l1/l2). See Keras Conv1D kernel_regularizer options for details.
    activity_regulariser  : str or None
        Activity regularisation to apply (e.g. l1, l2, l1/l2). See Keras Conv1D kernel_regularizer options for details.
    dropout               : float
        Dropout to apply before new CNN layer.

    Returns
    -------
    tf.keras.Model
        The Sequence UNET model
    """
    bottom_model = models.load_model(bottom_model, custom_objects=CUSTOM_OBJECTS)

    bottom_model.trainable = True

    if tune_layers > -1:
        for layer in bottom_model.layers[:-tune_layers]:
            layer.trainable = False

    input_features = bottom_model.inputs
    x = bottom_model.layers[-2 if features else -1].output

    if dropout:
        x = layers.SpatialDropout1D(dropout, name="top_dropout")(x)

    preds = layers.Conv1D(output_size, kernel_size, 1, padding='same',
                          activation=activation, name='preds',
                          kernel_regularizer=kernel_regulariser,
                          activity_regularizer=activity_regulariser)(x)

    return models.Model(inputs=input_features, outputs=preds)
