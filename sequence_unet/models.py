"""
Load, download and initialise trained and untrained Sequence UNET models.
"""
import requests
import os
from tensorflow.keras import layers, models

from sequence_unet.graph_cnn import GraphCNN
from sequence_unet import metrics

# TODO function to convert freq matrix to PSSM values?

__all__ = ["sequence_unet", "cnn_top_model", "download_trained_model",
           "download_all_models", "load_trained_model", "MODELS", "CUSTOM_OBJECTS"]

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

MODELS = {
	'freq_classifier': '',
	'pregraph_freq_classifier': '',
	'pssm_predictor': '',
	'pregraph_pssm_predictor': '',
	'patho_top': '',
	'pregraph_patho_top': '',
	'patho_finetune': '',
	'pregraph_patho_finetune': ''
}
"""
Dictionary mapping the IDs and download locations of each trained Sequence UNET model.
"""

def download_trained_model(model, root="."):
	"""
	Download a trained Sequence UNET model.

	Download the specified trained Sequence UNET model from BioStudies.
	Models are specified using the IDs indicated in models.MODELS, which also maps them to BioStudies files.
	Each model comes in a sequence only version (X) and version accepting additional structural input (pregraph_X).

	# Available models:

	* `freq_classifier`: Classifier predicting where variants occur above or below 0.01 observation frequency in a cross species multiple sequence alignment, as a proxy for deleteriousness.
	* `pssm_predictor`: Model predicting multiple alignment frequencies, which can be converted into a PSSM.
	* `patho_top`: Classifier predicting variant pathogenicity, trained as a new classifier head for `freq_classifier` on ClinVar data.
	* `patho_finetune`: Classifier predicting variant pathogenicity, trained by finetuning `freq_classifier` on ClinVar data.
	* `pregraph_freq_classifier`: Equivalent to `freq_classifier` taking structural input.
	* `pregraph_pssm_predictor`: Equivalent to `pssm_predictor` taking structural input.
	* `pregraph_patho_top`: Equivalent to `patho_top` taking structural input.
	* `pregraph_patho_finetune`: Equivalent to `patho_finetune` taking structural input.

	Parameters
	----------
	model : str
	        Sequence UNET model to download (see options in description and models.MODELS).
	root  : str
	        Root directory to download to.
	"""
	if not os.path.isdir(root):
		raise FileNotFoundError(f"No directory found at {root}")

	if model not in MODELS.keys():
		raise ValueError(f"model not recognised. Must one of: {', '.join(MODELS.keys())}")

	# download

def download_all_models(root="."):
	"""
	Download all trained Sequence UNET models.

	Download all trained Sequence UNET models from BioStudies.
	"""
	if not os.path.isdir(root):
		raise FileNotFoundError(f"No directory found at {root}")

	for model in MODELS.keys():
		download_trained_model(model, root=root)

def load_trained_model(model, root=".", download=False):
	"""
	Load a trained Sequence UNET model
	"""
	if os.path.exists(model):
		path = model

	elif model in MODELS.keys():
		path = f"{root}/{model}.tf"
		if not os.path.exists(path):
			if download:
				download_trained_model(model, root)
			else:
				raise FileNotFoundError("No model found at {path}")

	else:
		raise ValueError(("Unrecognised model - must be a path to a SavedModel or one of the "
		                  "options in sequence_unet.models.MODELS"))

	return models.load_model(path, custom_objects=CUSTOM_OBJECTS)

def sequence_unet(filters=8, kernel_size=5, num_layers=4, dropout=0,
                  graph_layers=None, graph_activation="relu",
                  conv_activation="relu", pred_activation="sigmoid",
                  kernel_regulariser=None, batch_normalisation=False):
    """
    Keras model predicting multiple alignment frequency matrix from sequences, based on
    a 1D UNET architecture.

    Expected input: (B, N, 20) arrays with B batches of Nx20 matrices, representing
                    sequences with one hot encoding on each row. An additional (B, N, N)
                    array with NxN contact graphs is required when graph_layers is not none.
                    The UNET architecture means N must be divisble by 2^(num_layers - 1).
                    This can be achieved with padding and masking.
    Output:         (B, N, 20) arrays of B batches of Nx20 arrays, with
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

    preds = layers.Conv1D(20, 1, 1, activation=pred_activation, name="predictor",
                          kernel_regularizer=kernel_regulariser)(x)

    return models.Model(inputs=inputs, outputs=preds)

def cnn_top_model(bottom_model=None, features=False, output_size=20,
                  tune_layers=3, kernel_size=3,
                  activation="sigmoid", kernel_regulariser=None,
                  activity_regulariser=None, dropout=0):
    """
    Add a 1D CNN layer top model to a trained Tensorflow Keras model

    bottom_model: PSSM model to train on top of. Bottom_model=None creates a simple CNN model.
    features: Start from model features rather than predicted output.
    tune_layers: Number of top layers to set as trainable.
    """
    if bottom_model is not None:
        bottom_model = models.load_model(bottom_model, custom_objects=CUSTOM_OBJECTS)

        bottom_model.trainable = True

        for layer in bottom_model.layers[:-tune_layers]:
            layer.trainable = False

        input_features = bottom_model.inputs
        x = bottom_model.layers[-2 if features else -1].output
    else:
        input_features = layers.Input(shape=[None, 20], name='input')
        x = input_features

    if dropout:
        x = layers.SpatialDropout1D(dropout, name="top_dropout")(x)

    preds = layers.Conv1D(output_size, kernel_size, 1, padding='same',
						  activation=activation, name='preds',
                          kernel_regularizer=kernel_regulariser,
                          activity_regularizer=activity_regulariser)(x)

    return models.Model(inputs=input_features, outputs=preds)
