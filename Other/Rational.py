class Rational:
    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __add__(self, r2):
        return self.add(r2)

    def __mul__(self, r2):
        return self.mult(r2)

    def invert(self):
        return Rational(-self.a, self.b)

    def add(self, r2):
        a1 = self.a
        b1 = self.b

        a2 = r2.a
        b2 = r2.b

        return Rational(a1*b2 + a2*b1, b1*b2)

    def mult(self, r2):
        a1 = self.a
        b1 = self.b

        a2 = r2.a
        b2 = r2.b

        return Rational(a1*b1, a2*b2)

    def __repr__(self):
        return f"{self.a}/{self.b}"
