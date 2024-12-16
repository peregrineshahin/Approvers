import re

file = open("..\src\search.c", "r")

pattern = r'^int\s+(\w+)\s*=\s*(\d+);'

for line in file.readlines():
    if re.search(pattern, line):
        (name, default) = re.search(pattern, line).groups()

        min = 1
        max = int(default) * 2
        step = round((min + max) / 20)
        if step == 0:
            step = 0.5

        print(f"{name}, int, {default}, {min}, {max}, {step}, 0.002")
