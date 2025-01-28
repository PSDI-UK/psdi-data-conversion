# convert.py
# Version 1.1, 31st August 2023
# Converts an SQLite database file to a byte array, and writes file Byte.js


data = []
array_contents = ''

try:
    with open("format.db", 'rb') as file:
        data = bytearray(file.read())
except Exception:
    print('An error occurred while opening or reading the file format.db')

for i in data:
    array_contents = array_contents + str(i) + ', '

array_contents = array_contents[:-2]

with open("byte.js", 'w') as file:
    file.write('const db_byte_array = [' + array_contents + '];' + '\n\n' + 'export default db_byte_array;')
