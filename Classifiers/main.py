import numpy as np
import sklearn.naive_bayes
from sklearn.neighbors import KNeighborsClassifier
import NBClassifier as nbc
import KnnClassifier as knn
import load_datasets as ld
import test_classifiers as tc


if __name__ == '__main__':
    dict = ld.load_knn_data_033_npy()
    X, y = dict['X_train'], dict['y_train']
    #X_test = dict['X_test']
    #Y_test = dict['Y_test']
    # X = np.array([[0, 0, 1],
    #               [0, 1, 1],
    #               [0, 0, 4],
    #               [0, 1, 4],
    #               [0, 0, 3],
    #               [0, 1, 3],
    #               [0, 0, 2],
    #               [0, 1, 2]], dtype=float)
    # y = np.array([1, 1, 4, 4, 4, 3, 2, 2])
    # X_test = np.array([[0, 0, 0], [0, 0, 2]], dtype=float)
    my_cls = knn.KnnBruteClassifier(weights='distance', n_neighbors=3)

    # my_cls.predict(X)
    # dist, _ = my_cls.kneighbors(X_test, 3)
    # tmp = my_cls._get_weights(dist)

    tc.test_classifier((nbc.MyGaussianNBClassifier(), sklearn.naive_bayes.GaussianNB(),
                         KNeighborsClassifier(), knn.KnnBruteClassifier(n_neighbors=5)), X, y)

