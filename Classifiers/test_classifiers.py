from sklearn.model_selection import train_test_split
import random


def test_classifier(classifiers, X, y, random_st=None, test_size=0.2):
    print("Compare {}-th classifiers".format(len(classifiers)))
    for i, classifier in enumerate(classifiers):
        print("classifier_{} {};".format(i, type(classifier).__name__))


    # Разделение данных на обучающую и тестовую выборки
    random_state = random.randrange(42) if random_st else 41
    print("random_state for train_test_split is ", random_state)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=random_state)

    print(f"train_size = {len(X_train)}\n"
          f"test_size  = {len(X_test)}")

    # Обучение классификатора
    for classifier in classifiers:
        classifier.fit(X_train, y_train)

    # Предсказание меток классов для тестовой выборки
    y_pred = []
    for classifier in classifiers:
        y_pred.append(classifier.predict(X_test))

    # Вывод результатов
    # for i, classifier in enumerate(classifiers):
    #     print("{} - labels predicted by {}".format(y_pred[i], type(classifier).__name__))
    # print("{} - true labels".format(y_test))
    for classifier in classifiers:
        print("Accuracy:", classifier.score(X_test, y_test))