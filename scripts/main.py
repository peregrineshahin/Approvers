import subprocess

p = "/kaggle_simulations/agent/approvers"
e = subprocess.Popen([p], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

def w(v):
    e.stdin.write(v)
    e.stdin.flush()

def main(v):
    t = v.remainingOverageTime * 1000
    w(f"stop\nkpos fen {v.board} moves {v.lastMove}\ngo wtime {str(t)} btime {str(t)}\n")
    return e.stdout.readline().strip().split()[1]
