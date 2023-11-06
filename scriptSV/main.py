import csv
import os
import numpy as np
import pandas as pd

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

        self.__stats_init(players_stats)
        self.__key_stats_init()
        self.__pass_stats_init()
        self.__def_stats_init()
        self.__strike_stats_init()

    def __stats_init(self, players_stats):
        dframe = pd.DataFrame(players_stats)
        stats_from_file = dframe.T
        stats_from_file.columns = stats_from_file.iloc[0]
        stats_from_file = stats_from_file.iloc[1:-1]
        stats_from_file = stats_from_file.reset_index(drop=True)
        stats_from_file.iloc[:, 1:] = stats_from_file.iloc[:, 1:].applymap(int)

        columns = (['Состав',
                    'Голы', 'Ассисты', 'Ключ. действ. : полож.', 'Ключ. действ.: отриц.',
                    'Короткий пас: точный', 'Короткий пас: неточный', 'Процент коротких',
                    'Длинный пас: точный', 'Длинный пас: неточный', 'Процент длинных', 'Выносы',
                    'Отбор: забрал', 'Отбор: В аут', 'Заработал на себе фолы', 'Сам сфолил',
                    'Удары в створ', 'Заблокированные', 'Мимо',
                    'Важное действие: полож.', 'Важное действие: отриц.', 'Гол со штрафного', 'Ключевое со штрафного',
                    'Не гол',
                    'Пенальти: гол', 'Пенальти: вратарь', 'Пенальти: мимо', 'Врт. пен: взял', 'Врт. пен: пропустил'])

        self.stats_df = pd.DataFrame(0, index=range(len(players_stats[0][1:-1])), columns=columns)
        number_of_teammates = self.stats_df.shape[0]
        self.stats_df['Состав'] = players_stats[0][1:-1]
        self.stats_df['Голы'] = stats_from_file['Goal - Sv'] if 'Goal - Sv' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Ассисты'] = stats_from_file['Асист'] if 'Асист' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Ключ. действ. : полож.'] = stats_from_file[
            'Ключевое действие - Положительное'] if 'Ключевое действие - Положительное' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Ключ. действ.: отриц.'] = stats_from_file[
            'Ключевое действие - Отрицательное'] if 'Ключевое действие - Отрицательное' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Короткий пас: точный'] = stats_from_file[
            'Пасы - Точный'] if 'Пасы - Точный' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Короткий пас: неточный'] = stats_from_file[
            'Пасы - Неточный'] if 'Пасы - Неточный' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Короткий пас: точный'] = pd.to_numeric(self.stats_df['Короткий пас: точный'], errors='coerce',
                                                              downcast='integer')
        self.stats_df['Короткий пас: неточный'] = pd.to_numeric(self.stats_df['Короткий пас: неточный'],
                                                                errors='coerce',
                                                                downcast='integer')
        self.stats_df['Процент коротких'] = (self.stats_df['Короткий пас: точный'] / (
                self.stats_df['Короткий пас: точный'] + self.stats_df['Короткий пас: неточный'])).mul(
            100).apply(lambda x: 0 if pd.isna(x) or x == float('inf') else x)
        self.stats_df['Длинный пас: точный'] = stats_from_file[
            'Длинные передачи - Точные'] if 'Длинные передачи - Точные' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Длинный пас: неточный'] = stats_from_file[
            'Длинные передачи - Неточные'] if 'Длинные передачи - Неточные' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Длинный пас: точный'] = pd.to_numeric(self.stats_df['Длинный пас: точный'], errors='coerce',
                                                             downcast='integer')
        self.stats_df['Длинный пас: неточный'] = pd.to_numeric(self.stats_df['Длинный пас: неточный'], errors='coerce',
                                                               downcast='integer')
        self.stats_df['Процент длинных'] = (self.stats_df['Длинный пас: точный'] / (
                self.stats_df['Длинный пас: точный'] + self.stats_df['Длинный пас: неточный'])).mul(
            100).apply(lambda x: 0 if pd.isna(x) or x == float('inf') else x)
        self.stats_df['Выносы'] = stats_from_file[
            'Длинные передачи - Вынос'] if 'Длинные передачи - Вынос' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Отбор: забрал'] = stats_from_file[
            'Отборы - Отобрал'] if 'Отборы - Отобрал' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Отбор: В аут'] = stats_from_file[
            'Отборы - В аут'] if 'Отборы - В аут' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Заработал на себе фолы'] = stats_from_file[
            'Фолы - Заработал'] if 'Фолы - Заработал' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Сам сфолил'] = stats_from_file[
            'Фолы - Сфолил'] if 'Фолы - Сфолил' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Удары в створ'] = stats_from_file[
            'Удары - В створ'] if 'Удары - В створ' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Заблокированные'] = stats_from_file[
            'Удары - Блок'] if 'Удары - Блок' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)
        self.stats_df['Мимо'] = stats_from_file['Удары - Мимо'] if 'Удары - Мимо' in stats_from_file else pd.Series(
            np.zeros(number_of_teammates), dtype=int)

        # need to do
        # Важное действие: полож.', 'Важное действие: отриц.', 'Гол со штрафного', 'Ключевое со штрафного',
        # 'Не гол',
        # 'Пенальти: гол', 'Пенальти: вратарь', 'Пенальти: мимо', 'Врт. пен: взял', 'Врт. пен: пропустил'])

        # print(self.stats_df.columns)
        # new_names = {'': '',
        #              'Goal - Sv': '',
        #              'Асист': '',
        #              'Substitution': '',
        #              'Ключевое действие - Положительное': '',
        #              'Ключевое действие - Отрицательное': '',
        #              'Отборы - Отобрал': '',
        #              'Отборы - В аут': '',
        #              'Пасы - Точный': '',
        #              'Пасы - Неточный': '',
        #              'Длинные передачи - Точные': '',
        #              'Длинные передачи - Вынос': '',
        #              'Длинные передачи - Неточные': '',
        #              'Удары - В створ': '',
        #              'Удары - Блок': '',
        #              'Удары - Мимо': '',
        #              'Фолы - Заработал': '',
        #              'Фолы - Сфолил': '',
        #              'Штрафные - Мимо': '',
        #              'Вратарь - Сейв': '',
        #              'Дриблинг - Успешно': '',
        #              'Дриблинг - Потеря': '',
        #              'TOTAL': 'TOTAL',
        #              'Mins Played': ''}

    def __key_stats_init(self):
        self.key_stats = self.stats_df[['Состав', 'Голы', 'Ассисты', 'Ключ. действ. : полож.', 'Ключ. действ.: отриц.']]
        self.key_stats = self.key_stats.sort_values(
            by=['Голы', 'Ассисты', 'Ключ. действ. : полож.', 'Ключ. действ.: отриц.'], ascending=False)

    def __pass_stats_init(self):
        self.pass_stats = self.stats_df[['Состав', 'Короткий пас: точный', 'Короткий пас: неточный', 'Процент коротких',
                                         'Длинный пас: точный', 'Длинный пас: неточный', 'Процент длинных', 'Выносы']]
        self.pass_stats.loc[:, 'Процент коротких'] = self.pass_stats['Процент коротких'].round(2)
        self.pass_stats.loc[:, 'Процент длинных'] = self.pass_stats['Процент длинных'].round(2)
        self.pass_stats = self.pass_stats.sort_values(
            by='Процент коротких', ascending=False)

    def __def_stats_init(self):
        self.def_stats = self.stats_df[
            ['Состав', 'Отбор: забрал', 'Отбор: В аут', 'Заработал на себе фолы', 'Сам сфолил']]
        self.def_stats = self.def_stats.sort_values(by='Состав')

    def __strike_stats_init(self):
        self.strike_stats = self.stats_df[['Состав', 'Удары в створ', 'Заблокированные', 'Мимо', 'Голы']]
        self.strike_stats = self.strike_stats.sort_values(by='Состав')

    @staticmethod
    def create_pdf_table(data, pdf_filename):
        # Создайте новое изображение Matplotlib
        fig, ax = plt.subplots()  # Измените размер в соответствии с вашими потребностями

        # Создайте таблицу с данными

        table = ax.table(cellText=data.values, colLabels=data.columns, loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(8)
        table.auto_set_column_width(range(data.shape[1]))

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


if __name__ == '__main__':
    with open('statistic.csv', 'r', encoding='utf8') as csvfile:
        reader = list(csv.reader(csvfile, delimiter=';'))
        players_statistic = reader[reader.index(['Actions/Players']) + 1:]
        print(players_statistic)
        stats = MatchAnalyze(players_statistic)
        print('Done')

        MatchAnalyze.create_pdf_table(stats.key_stats, "key_stats.pdf")
        MatchAnalyze.create_pdf_table(stats.pass_stats, "pass_stats.pdf")
        MatchAnalyze.create_pdf_table(stats.def_stats, "def_stats.pdf")
        MatchAnalyze.create_pdf_table(stats.strike_stats, "strike_stata.pdf")
