import numpy as np


class KnnBruteClassifier(object):
    """Классификатор реализует взвешенное голосование по ближайшим соседям.
    Поиск ближайшего соседа осуществляется полным перебором.
    Параметры
    ----------
    n_neighbors : int, optional
        Число ближайших соседей, учитывающихся в голосовании
    weights : str, optional (default = 'uniform')
        веса, используемые в голосовании. Возможные значения:
        - 'uniform' : все веса равны.
        - 'distance' : веса обратно пропорциональны расстоянию до классифицируемого объекта
        - функция, которая получает на вход массив расстояний и возвращает массив весов
    metric: функция подсчета расстояния (по умолчанию l2).
    """

    def __init__(self, n_neighbors=1, weights='uniform', metric="l2"):
        self.classes = None
        self.X = None
        self.Y = None
        self.n_neighbors = n_neighbors
        self.weights = weights
        self.metric = metric

    def fit(self, x, y):
        """Обучение модели.
        Параметры
        ----------
        x : двумерным массив признаков размера n_queries x n_features
        y : массив/список правильных меток размера n_queries
        Выход
        -------
        Метод возвращает обученную модель
        """
        if x.ndim != 2 or y.ndim != 1 or len(x) != len(y):
            raise ValueError('Wrong data X or y')
        self.X = np.copy(x)
        self.Y = np.copy(y)
        self.classes = np.unique(self.Y)
        return self

    def predict(self, X):
        """ Предсказание класса для входных объектов
        Параметры
        ----------
        X : двумерным массив признаков размера n_queries x n_features
        Выход
        -------
        y : Массив размера n_queries
        """
        return self.classes[np.argmax(self.predict_proba(X), axis=1)]

    def predict_proba(self, X):
        """Возвращает вероятности классов для входных объектов
        Параметры
        ----------
        X : двумерным массив признаков размера n_queries x n_features
        Выход
        -------
        p : массив размера n_queries x n_classes] c вероятностями принадлежности
        объекта к каждому классу
        """
        predicts = np.zeros((len(X), len(self.classes)))
        neigh_dist, neigh_ind = self.kneighbors(X, self.n_neighbors)
        neigh_weights = self._get_weights(neigh_dist)
        neigh_classes = self.Y[neigh_ind]

        for i, predict in enumerate(predicts):
            for k, class_k in enumerate(self.classes):
                predict[k] = np.sum(neigh_weights[i][np.where(neigh_classes[i] == class_k)])

        return predicts

    def kneighbors(self, X, n_neighbors):
        """Возвращает n_neighbors ближайших соседей для всех входных объектов и расстояния до них
        Параметры
        ----------
        X : двумерным массив признаков размера n_queries x n_features
        Выход
        -------
        neigh_dist массив размера n_queries х n_neighbors
        расстояния до ближайших элементов
        neigh_indarray, массив размера n_queries x n_neighbors
        индексы ближайших элементов
        """
        neigh_dist = np.zeros((len(X), n_neighbors), dtype=float)
        neigh_indarray = np.zeros((len(X), n_neighbors), dtype=int)

        for i, x in enumerate(X):
            closes_neigh = np.argpartition(np.linalg.norm(self.X - x, axis=1), n_neighbors)[:n_neighbors]
            neigh_indarray[i] = closes_neigh
            neigh_dist[i] = np.linalg.norm(self.X[closes_neigh] - x, axis=1)

        return neigh_dist, neigh_indarray

    def _get_weights(self, dist):
        if callable(self.weights):
            return self.weights(dist)
        if self.weights == 'distance':
            return np.divide(1, dist, out=np.full(dist.shape, np.inf), where=dist != 0)
        if self.weights == 'uniform':
            return np.ones(dist.shape)
        raise ValueError("Weights should be one from: 'uniform', 'distance' or callable'")

    def score(self, X, y):
        return np.mean(self.predict(X) == y)
