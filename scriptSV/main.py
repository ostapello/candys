import csv
import os
import numpy as np
import pandas as df

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tabulate import tabulate
import pdfkit


class MatchAnalyze:
    def __init__(self, players_stats):
        self.all_table = players_stats
        self.team = players_stats[0][1:-1]
        self.stats_name = [stat[0] for stat in players_stats[1:]]
        self.stats = np.zeros((len(self.stats_name), len(self.team)), dtype='int8')
        for i, line in enumerate(players_stats[1:]):
            self.stats[i] = np.array(line[1:-1])

        self.__key_stats_init()

    def __key_stats_init(self):
        # Формирую статистику для ключевых действий: Голы, Асисты, Ключевые действия
        key_stats = ['Goal - Sv', 'Асист', 'Ключевое действие - Положительное', 'Ключевое действие - Отрицательное']
        table_key_stats = [['Состав', 'Голы', 'Ассисты', 'Полож. ключ. действ.', 'Отриц. ключ. действ.']]
        for i, name in enumerate(self.team):
            row = [name]
            for stat in key_stats:
                try:
                    index = self.stats_name.index(stat)
                except:
                    index = None
                if index is not None:
                    row.append(str(self.stats[index][i]))
                else:
                    row.append('-')
            table_key_stats.append(row)
        self.key_stats = df.DataFrame(data=table_key_stats[1:], columns=table_key_stats[0])
        self.key_stats = self.key_stats.sort_values(by='Голы', ascending=False)

    @staticmethod
    def create_pdf_table(data, pdf_filename):
        # Создайте новое изображение Matplotlib
        fig, ax = plt.subplots()  # Измените размер в соответствии с вашими потребностями

        # Создайте таблицу с данными
        table_data = []
        for row in data:
            table_data.append(row)

        table = ax.table(cellText=table_data, loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.auto_set_column_width(range(len(data[0])))

        ax.axis('off')  # Скрыть оси координат

        # Показать рамки таблицы
        for key, cell in table._cells.items():
            cell.set_linewidth(0.5)  # Установите толщину рамки в соответствии с вашими потребностями

        # Создайте новый PDF-файл
        pdf_pages = PdfPages(pdf_filename)

        # Добавьте таблицу на страницу PDF
        pdf_pages.savefig(fig, bbox_inches='tight')

        # Закройте PDF-файл
        pdf_pages.close()

    @staticmethod
    def save_df_as_pdf(dataframe, pdf_filename):
        # Используйте tabulate для форматирования DataFrame
        table = tabulate(dataframe, headers='keys', tablefmt='fancy_grid')

        # Создайте временный HTML-файл с таблицей
        html_filename = 'temp_table.html'
        with open(html_filename, 'w', encoding='utf8') as html_file:
            html_file.write(table)

        # Преобразуйте HTML-файл в PDF с помощью pdfkit
        pdfkit.from_file(html_filename, pdf_filename)

        # Удалите временный HTML-файл
        os.remove(html_filename)


if __name__ == '__main__':
    with open('statistic.csv', 'r', encoding='utf8') as csvfile:
        reader = list(csv.reader(csvfile, delimiter=';'))
        players_statistic = reader[reader.index(['Actions/Players']) + 1:]
        print(players_statistic)
        stats = MatchAnalyze(players_statistic)
        print('Done')
        # Пример использования:
        MatchAnalyze.create_pdf_table(stats.key_stats, "output.pdf")
        MatchAnalyze.save_df_as_pdf(stats.key_stats, "output1.pdf")
