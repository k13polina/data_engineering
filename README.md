# data_engineering
## Это репозиторий под анализ дифференциальной экспрессии генов на модельном растении _Arabidopsis thaliana_
### На входе у нас есть данные с прибора, в который поместили образцы РНК, выделенной из листьев растений, собранных при разных условиях, в ИППИ РАН.
### 1. Результаты проверки качества первичных данных находятся в папке **reports_fastqc** и свидетельствуют о хорошем качестве.
### 2. Для предобработки данных написан пайплайн - __script_preprocessing.bash__. Полученные результаты внесены в базу данных __genes_activity.db__.
### 3. Для дальнейшего исследования написан пайплайн на языке R - __script_processing.R__. Результаты в папке **reports_prosessing**.
### 4. Провели предварительный анализ функций генов с помощью Gene Ontology. 
### [Дашборд](https://datalens.yandex/xorku7gnnwuaj) по полученным данным представлен по ссылке.
#### Выявили, что наиболее часто встречаются гены, дифференциально экспрессирующиеся в ходе циркадного цикла, связанного с фотосинтезом, так как изменяющимся фактором является количество солнечного света, которое влияет на продуктивность данного цикла.
