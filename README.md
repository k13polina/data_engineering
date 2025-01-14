# data_engineering
## Это репозиторий итогового проекта по курсу «Инжинеринг управления данными» Лиходзиевской Маргариты и Кузнецовой Полины на тему «Повышение продуктивности сельскохозяйственных растений за счёт анализа транскриптомных данных».
## Презентация с подробным описанием хода работы и визуальзацией результатов по [ссылке](https://docs.google.com/presentation/d/1WZt6jD-K0fF-tKCVHfnt9IXpx8MOaQhQlHbGA_BlH08/edit#slide=id.p)
На входе у нас есть данные с прибора, в который поместили образцы РНК, выделенной из листьев модельного растения _Arabidopsis thaliana_, собранных в ИППИ РАН при разных условиях.
1. Результаты проверки качества первичных данных находятся в папке [**reports_fastqc**](https://github.com/k13polina/data_engineering/tree/main/reports_fastqc) и свидетельствуют о хорошем качестве.
2. Для предобработки данных написан пайплайн - [__script_preprocessing.bash__](https://github.com/k13polina/data_engineering/blob/main/script_preprocessing.bash). Полученные результаты внесены в базу данных __genes_activity.db__.
3. Для дальнейшего исследования написан пайплайн на языке R - [__script_processing.R__](https://github.com/k13polina/data_engineering/blob/main/script_processing.R). Результаты в папке [**reports_prosessing**](https://github.com/k13polina/data_engineering/tree/main/reports_prosessing).
4. Провели предварительный анализ функций генов с помощью Gene Ontology. 
Дашборд по полученным данным представлен по [ссылке](https://datalens.yandex/xorku7gnnwuaj).
Выявили, что наиболее часто встречаются гены, дифференциально экспрессирующиеся в ходе циркадного цикла, связанного с фотосинтезом, так как изменяющимся фактором является количество солнечного света, которое влияет на продуктивность данного цикла.
