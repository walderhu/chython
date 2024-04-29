from chython import smiles
import time
import sys
database_smiles = 'smiles_data.txt'

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

with open(database_smiles, 'r') as smiles_file, open('errors.txt', 'a') as errors_file:
    lines = smiles_file.readlines()
    Runs = len(lines)
    for i, line in enumerate(lines):
        smile = line
        start_time = time.time() 
        try:
            m = smiles(smile)
            m.clean2d()
        except Exception as e:
            errors_file.write(f"{smile} - {str(e)}\n")
        else:
            elapsed_time = time.time() - start_time
            if elapsed_time > 5:
                errors_file.write(f"{smile} - Processing took more than 5 seconds\n")
        
        progressBar = "\rProgress: " + ProgressBar(Runs, i + 1, BarLength=20, ProgressIcon="#", BarIcon="-")
        ShowBar(progressBar)
        
print('\nDone.')
