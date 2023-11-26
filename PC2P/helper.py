from termcolor import colored

def printc(str, *args, **kwargs):
    print(colored(str, "red"), end = kwargs.get('end', None))