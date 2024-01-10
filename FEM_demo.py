import subprocess
folder_path = r"C:\Users\theco\PycharmProjects\ZirconiaStage"
script_path = folder_path+"\poisson_neumann.py"
subprocess.run(['cd', '/D', folder_path], shell=True)
subprocess.run(['sfepy-run', script_path], cwd=folder_path, shell=True)

