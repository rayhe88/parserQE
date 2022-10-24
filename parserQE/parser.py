
def loadFile(name):
    with open(name) as file:
        #for key, value in iter(file.readlines()[::-1]):
        for key, value in enumerate(iter(file.readlines()[::-1])):
            print(key, value.strip())
            if key == 10:
                break
