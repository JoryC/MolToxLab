import os


def clean_file(filename):
  with open(filename, 'r') as f:
    s = f.read()
    s = s[s.find('Probe'):]

  with open(filename, 'w') as f:
    f.write(s)
 

def getfiles():
  listing = os.listdir()
  txts = []

  for file in listing:
    if file.endswith('.txt'):
      txts.append(file)

  return txts


def main():
  files = getfiles()

  for file in files:
    clean_file(file)


if __name__ == '__main__':
  main()
