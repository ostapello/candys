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
        -  функция, которая получает на вход массив расстояний и возвращает массив весов
    metric: функция подсчета расстояния (по умолчанию l2).
    """

    def __init__(self, n_neighbors=1, weights='uniform', metric="l2"):
        pass

    def fit(self, x, y):
        """Обучение модели.
        Парметры
        ----------
        x : двумерным массив признаков размера n_queries x n_features
        y : массив/список правильных меток размера n_queries
        Выход
        -------
        Метод возвращает обученную модель
        """
        pass

    def predict(self, x):
        """ Предсказание класса для входных объектов
        Параметры
        ----------
        X : двумерным массив признаков размера n_queries x n_features
        Выход
        -------
        y : Массив размера n_queries
        """
        pass

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
        pass

    def kneighbors(self, x, n_neighbors):
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
        pass