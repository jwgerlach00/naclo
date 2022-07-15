import argparse


argparser = argparse.ArgumentParser(description='Convert a file from one encoding to another.')
argparser.add_argument('-f', '--file', help='Input file path.', required=True)
path = argparser.parse_args().file

with open(path) as f:
    data = f.read()

data = data.replace('"', '\'').replace('true', 'True').replace('false', 'False')

print(data)