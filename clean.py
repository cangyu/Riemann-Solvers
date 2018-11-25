import os

def action():
    os.system('rm *.txt')
    os.system('rm *.out')


def clean(root_dir):
    t = os.path.abspath(root_dir)
    t = os.path.split(t)[-1]
    if(t == '.git'):
        return
    
    action()
    for f in os.listdir(root_dir):
        if os.path.isdir(f):
            os.chdir(os.path.join('.', f))
            clean('.')
            os.chdir('..')

if __name__ == '__main__':
    clean(os.getcwd())