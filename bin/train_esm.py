#!/usr/bin/env python3
"""
Train PSSM and frequency classification top models based on ESM1b using sequences in a ProteinNet file.
"""
import sys
import argparse
import torch
import proteinnetpy
import numpy as np

class ESM1bData(torch.utils.data.IterableDataset):
    """Iterable dataset returning ESM1b predictions from a ProteinNet dataset"""
    def __init__(self, proteinnet, no_gpu=False, thresh=None):
        super(ESM1bData).__init__()
        self.thresh = thresh
        self.use_gpu = torch.cuda.is_available() and not no_gpu

        # Load ESM1b modek
        self.model, self.alphabet = torch.hub.load("facebookresearch/esm:main",
                                                   "esm1b_t33_650M_UR50S")
        self.batch_converter = self.alphabet.get_batch_converter()
        self.model = self.model.eval()

        if self.use_gpu:
            self.model = self.model.cuda()
            print("Transferred ESM1b to GPU", file=sys.stderr)

        # Load ProtienNet data
        print(f"Initialising ProtienNet dataset: {proteinnet}", file=sys.stderr)
        filter_func = proteinnetpy.data.make_length_filter(max_length=1022)
        self.proteinnet = proteinnetpy.data.ProteinNetDataset(path=proteinnet,
                                                              filter_func=filter_func,
                                                              preload=False)

        # Length calcualted on init to not have to re-stream through ProteinNet file
        self.length = sum(len(r) for r in self.proteinnet)

    def __len__(self):
        return self.length

    def __iter__(self):
        return self.generate_embeddings()

    def generate_embeddings(self):
        """
        Generator yielding per position embeddingsa and corresponding PSSMs or frequency categories
        """
        for record in self.proteinnet:
            esm_input = [(record.id, "".join(record.primary))]
            _, _, batch_tokens = self.batch_converter(esm_input)

            if self.use_gpu:
                batch_tokens = batch_tokens.to(device="cuda", non_blocking=True)

            with torch.no_grad():
                results = self.model(batch_tokens, repr_layers=[33], return_contacts=True)
            reps = results["representations"][33].to(device="cpu")

            for p in range(len(record)):
                x = reps.numpy()[p,:]
                y = record.evolutionary[:,p]

                if self.thesh is not None:
                    y = (y < self.thresh).astype(float)

                yield x, y

class ESM1bMock(torch.utils.data.IterableDataset):
    """Iterable dataset returning ESM1b predictions from a ProteinNet dataset"""
    def __init__(self, thresh=None):
        super(ESM1bMock).__init__()
        self.thresh = thresh

    def __len__(self):
        return 100

    def __iter__(self):
        return self.generate_embeddings()

    def generate_embeddings(self):
        """
        Generator yielding per position embeddingsa and corresponding PSSMs or frequency categories
        """
        for _ in range(100):
            x = np.random.normal(size=1280)
            y = np.random.random(size=20) if self.thresh is None else np.random.randint(2, size=20)
            yield x, y

class LinearModel(torch.nn.Module):
    """Simple linear top layer"""
    def __init__(self, classifier=False):
        super(LinearModel, self).__init__()
        self.linear = torch.nn.Linear(1280, 20)
        self.activation = torch.nn.Softmax(dim=1) if classifier else torch.nn.Sigmoid()

    def forward(self, x):
        out = self.linear(x)
        out = self.activation(out)
        return out

def main():
    """
    Run ESM1b on input Fasta and format into a single output file
    """
    args = parse_args()
    classifier = args.threshold is not None

    dataloader = torch.utils.data.DataLoader(ESM1bData(proteinnet=args.proteinnet,
                                                       no_gpu=args.nogpu,
                                                       thresh=args.threshold),
                                             batch_size=1000)

    model = LinearModel(classifier=classifier)

    loss_fn = torch.nn.BCELoss() if classifier else torch.nn.MSELoss()
    optimiser = torch.optim.SGD(model.parameters(), lr=1e-3)

    size = len(dataloader.dataset)
    model.train()
    for epoch in range(5):
        print(f"Epoch {epoch+1}\n-------------------------------")
        for batch, (x, y) in enumerate(dataloader):
            x, y = x.float(), y.float()

            # Compute prediction error
            pred = model(x)
            loss = loss_fn(pred, y)
            
            # Backpropagation
            optimiser.zero_grad()
            loss.backward()
            optimiser.step()
            if batch % 100 == 0:
                loss, current = loss.item(), batch * len(x)
                print(f"loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")

    torch.save(model, args.model)




def parse_args():
    """Process arguments"""
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('proteinnet', metavar='P', type=str, help="Input ProteinNet data")

    parser.add_argument('--no_gpu', '-n', action="store_true",
                        help="Prevent GPU usage for ESM1b even when available")

    parser.add_argument('--threshold', '-r', default=None, type=float,
                        help="Perform frequency classification at given threshold")

    parser.add_argument('--model', '-m', default="esm_top_model.pth", help="Path to save model")

    parser.add_argument('--epochs', '-e', default=3, type=int, help="Epochs to train for")

    return parser.parse_args()

if __name__ == "__main__":
    main()
