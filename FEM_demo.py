file_path = 'OutputLog.txt'
with open(file_path, 'r') as file:
    lines = file.readlines()
lines = [line for line in lines if line.startswith('Comput')]
with open(file_path, 'w') as file:
    file.writelines(lines)
