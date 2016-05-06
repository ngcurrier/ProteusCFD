import random

class colors():
    def __init__(self):
        self.red = (1,0,0)
        self.blue = (0,0,1)
        self.yellow = (1,1,0)
        self.aqua = (0,1,1)
        self.magenta = (1,0,1)
        self.gray = (0.5, 0.5, 0.5)
        self.white = (0,0,0)
        self.black = (1,1,1)
        random.seed(999.0)
        self.list = []
        self.list.append(self.red)
        self.list.append(self.blue)
        self.list.append(self.yellow)
        self.list.append(self.aqua)
        self.list.append(self.magenta)
        self.list.append(self.white)
        self.list.append(self.gray)
        self.list.append(self.black)
        print self.list
        self.counter = 0
        
    def getRandom(self):
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()
        return (r1,r2,r3)

    def getNext(self):
        if self.counter > len(self.list):
            raise 'colors.getNext() list too short'
        color = self.list[self.counter]
        self.counter = self.counter + 1
        return color
