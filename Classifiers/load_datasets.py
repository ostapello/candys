import numpy as np
from pathlib import Path
from sklearn.datasets import load_iris


def load_mnist_npz(path=Path('datasets').joinpath('mnist.npz')):
    print(path)
    with np.load(path, allow_pickle=True) as f:  # pylint: disable=unexpected-keyword-arg
        x_train, y_train = f['x_train'], f['y_train']
        x_test, y_test = f['x_test'], f['y_test']

    return (x_train, y_train), (x_test, y_test)


def load_knn_data_033_npy(path=Path('datasets').joinpath('knn_data_033.npy')):
    return np.load(path, allow_pickle=True).tolist()


def load_skl_iris():
    return load_iris(return_X_y=True)
