import pandas as pd

# Чтение данных из файла Excel
df = pd.read_excel('Dataset_and_Predictions.xlsx')

# Извлечение данных из столбца 'SMILES'
smiles_data = df['SMILES']

# Сохранение данных в текстовый файл
with open('smiles_data.txt', 'w') as file:
    for smiles in smiles_data:
        file.write(f"{smiles}\n")