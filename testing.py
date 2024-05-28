from chython import smiles
import time
import sys
import os
# import shutil
def main():
    database_smiles = 'smiles_database.txt'
    normal_errors_file = 'errors.txt'
    time_errors_file = 'time_error.txt'
    output_folder = 'output_files'  # Папка для вывода файлов

    # Создаем папку для вывода, если она еще не существует
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    def ProgressBar(Total, Progress, BarLength=20, ProgressIcon="#", BarIcon="-"):
        try:
            if BarLength < 1:
                BarLength = 20
            Status = ""
            Progress = float(Progress) / float(Total)
            if Progress >= 1.:
                Progress = 1
                Status = "\r\n"
            Block = int(round(BarLength * Progress))
            Bar = "[{}] {:.1f}% {}".format(ProgressIcon * Block + BarIcon * (BarLength - Block), Progress * 100, Status)
            return Bar
        except:
            return "ERROR"

    def ShowBar(Bar):
        sys.stdout.write(Bar)
        sys.stdout.flush()

    with open(database_smiles, 'r') as smiles_file, \
        open(os.path.join(output_folder, normal_errors_file), 'a') as normal_errors_file, \
        open(os.path.join(output_folder, time_errors_file), 'a') as time_errors_file:
        lines = smiles_file.readlines()
        Runs = len(lines)
        for i, line in enumerate(lines):
            smile = line.strip('\n')
            start_time = time.time()
            try:
                m = smiles(smile)
                m.clean2d()
            except Exception as e:
                error_info = f"{smile}\n"
                normal_errors_file.write(error_info)
            else:
                elapsed_time = time.time() - start_time
                if elapsed_time > 5:
                    error_info = f"{smile}\n"
                    time_errors_file.write(error_info)
            
            progressBar = "\rProgress: " + ProgressBar(Runs, i + 1, BarLength=20, ProgressIcon="#", BarIcon="-")
            ShowBar(progressBar)

    print('\nDone.')

if __name__ == '__main__':
    main()