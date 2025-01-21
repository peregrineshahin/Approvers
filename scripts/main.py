import subprocess

p = "/kaggle_simulations/agent/cfish"
e = subprocess.Popen([p], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

def w(v):
    e.stdin.write(v)
    e.stdin.flush()

def main(v):
    t = v.remainingOverageTime * 1000
    w(f"stop\nposition fen {v.board}\ngo wtime {str(t)} btime {str(t)}\n")
    return e.stdout.readline().strip().split()[1]
