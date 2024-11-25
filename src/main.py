import os
import subprocess

class Engine:

    def __init__(self, engine_file):
        self.engine_file = engine_file
        self.bestmove = ""
        self.engine_process = None

    def write(self, input_text):
        self.engine_process.stdin.write(input_text)
        self.engine_process.stdin.flush()

    def connect(self):
        self.engine_process = subprocess.Popen([self.engine_file], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                               universal_newlines=True)

        self.write("uci\n")
        self.write("ucinewgame\n")
        self.write("isready\n")

    def listen(self):
        while True:
            response = self.engine_process.stdout.readline().strip()
            tokens = response.split()

            if len(tokens) >= 2 and tokens[0] == "bestmove":
                self.bestmove = tokens[1]
                break

    def think(self, allocated_time, fen):
        command = "position fen " + fen + "\n"

        self.write(command)
        self.write("isready\n")
        self.write("go wtime " + str(allocated_time) + " btime " + str(allocated_time) + "\n")

        self.listen()


engine_name = "cfish_compressed"
engine_file_path = "./" + engine_name

if not os.path.exists(engine_file_path):
    engine_file_path = "/kaggle_simulations/agent/" + engine_name
elif not os.path.exists(engine_file_path):
    engine_file_path = "/kaggle/working/" + engine_name
elif not os.path.exists(engine_file_path):
    engine_file_path = "/kaggle/working/src/" + engine_name

engine = Engine(engine_file_path)
engine.connect()


def main(obs):
    fen = obs.board
    if "r1bQk2r/4ppbp/p1p2np1/4p3/2P5/2N5/PP2BPPP/R1B2RK1" in fen:
        return "e8d8"
    time_left = obs.remainingOverageTime * 1000

    engine.think(time_left, fen)

    return engine.bestmove
