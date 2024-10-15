import math

# Problem 1

def problem1():
    """Prompt user for height, then calculate time it takes to fall"""
    print("--- PROBLEM 1 ---\n")
    print("A ball is falling from a tower.")
    h = int(input("How high is the tower in meters? "))

    g = 9.80 # m/s^2 

    time = math.sqrt((2 * h) / g)

    print(f"It would take {time = } seconds to fall a distance of {h} meters.")

def problem2():
    """Print all the Catalan numbers under 1 billion"""
    print("\n--- PROBLEM 2 ---\n")

    largest = 1_000_000_000

    n = 0 # store the index of the current Catalan number
    current = 1 # store the current Catalan number

    while current < largest:
        print(current)
        current = int((4 * n + 2)/(n + 2) * current)
        n += 1

def problem3():
    """Print the first 20 lines of Pascal's Triangle"""
    print("\n--- PROBLEM 3 ---\n")

    # part (a)
    def factorial(n: int) -> int:
        """Returns n!"""
        if n == 0:
            return 1
        else:
            return n * factorial(n-1)

    def binomial(n: int, k: int) -> int:
        """Returns n choose k"""
        return int(factorial(n) / (factorial(k) * factorial(n - k)))

    # part (b)
    current_row = ""
    for n in range(1, 21):
        for k in range(0, n+1):
            current_row += str(binomial(n, k))
        print(current_row)
        current_row = ""

def problem4():
    """Recursion problems."""
    print("\n--- PROBLEM 4 ---\n")

    # part (a)
    def catalan(n: int) -> int:
        """Returns the nth Catalan number, according to the given formula"""
        if n == 0:
            return 1
        else:
            # the given formula was for n+1, so I just subtract 1 to get for n
            return int((4*(n-1) + 2)/((n-1) + 2) * catalan(n-1))

    print(f"The 100th Catalan number is {catalan(100)}")

    # part (b)
    def g(m: int, n: int) -> int:
        """Returns the greatest common divisor of m and n."""
        if n == 0:
            return m
        else:
            return g(n, m % n)

    print(f"The greatest common divisor of 108 and 192 is {g(108, 192)}")

def problem5():
    """Prime numbers"""
    print("\n--- PROBLEM 5 ---\n")
    limit = 10_000
    primes = [2]

    for i in range(3, limit+1):
        is_prime = True # boolean to store the current status of i

        for prime in primes: # look at each existing prime
            if prime >= math.sqrt(i):
                break # stop looking at primes 
            if i % prime == 0:
                is_prime = False # i is divisible by this prime, so i is not prime

        if is_prime:
            primes.append(i)

    print(f"The primes under {limit} are: {primes}")

if __name__ == "__main__":
    problem1()
    problem2()
    problem3()
    problem4()
    problem5()
