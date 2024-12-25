import os
import subprocess

p = "/kaggle_simulations/agent/cfish"
os.system("chmod +x " + p)
e = subprocess.Popen([p], stdin=subprocess.PIPE, stdout=subprocess.PIPE, universal_newlines=True)

def w(v):
    e.stdin.write(v)
    e.stdin.flush()

def main(v):
    t = v.remainingOverageTime * 1000
    w(f"stop\nposition fen {v.board}\ngo time {str(t)}\n")
    return e.stdout.readline().strip().split()[1]
