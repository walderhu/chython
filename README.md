
### chython

![all_err](misc/images/all_err.png)

Пропустил через проду файл smiles_database.txt с более чем 12к примеров smiles, 5789 завершили работу с ошибкой (подробно можно ознакомиться в файле all_errors.txt), категории ошибок в основном делятся на 3 типа:
- не реализован алгоритм обработки явного отображения водородных связей (самая распространенная ошибка)
- возникает ошибка кекулизации молекулы
- ошибка типа (встречается крайне редко)
2 и 3 ошибки возникают скорее всего тоже изза неправильной обработки явного отображения водородных связей

В файле main.ipynb можно посмотреть основные ошибки и примеры работы действующей программы

в файле testing.py проводиться тестировка и прогонка через программу всех примеров (всех 12 тысяч), довольно долгий в сумме процесс, если хотите его проверить самостоятельно, его можно запустить через python testing.py или более удобным вызовом через makefile путем ввода в терминал цели make test

в файле unit реализованы функциональные юнит тесты на ошибки обработки программой smile-строки. Самостоятельно запустить юнит тесты и посмотреть репорт в виде html(на используемость функций написанных мною на реальных примерах) файла можно через команду make gcov. Вариант моего запуска приложен ниже.
![htmlcov](misc/images/htmlcov.png)

В настоящее время ведется реализация и оптимизация алгоритма обработки отрисовки явных водородов с проекта SmileDrawwer, который написанный на js. (проблемма в том, что этот алгоритм явно не вынесен в отдельные функции, а разбросан по всему проекту, поэтому нужно проводить аналитику отработки кода на js, дебаггера или gcov тестов для js я не знаю, поэтому приходится отслеживать все самому).

Также после реализации и тестирования этого алгоритма, нужно будет расписать проведенный форк проекта пикачу, рассписать всю проведенную работу (по просьбе Валентины Александровны).