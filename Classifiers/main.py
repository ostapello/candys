import sklearn.naive_bayes
import sklearn.neighbors
import NBClassifier as nbc
import KnnClassifier as knn
import load_datasets as ld
import test_classifiers as tc


if __name__ == '__main__':
    dict = ld.load_knn_data_033_npy()
    X, y = dict['X_train'], dict['y_train']
    tc.test_classifier((nbc.MyGaussianNBClassifier(), sklearn.naive_bayes.GaussianNB(),
                        sklearn.neighbors.KNeighborsClassifier()), X, y)

