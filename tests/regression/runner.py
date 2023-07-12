# Script to run regression tests and check output

import subprocess

# Make Emplode executable
subprocess.run(["make", "quick"], cwd="../..")
print()

success = 0
failure = 0

for test in ["hello_world"]:
    # Find expected output
    file = open(test + ".emp", "r")
    line = file.readline()
    assert(line.startswith("// Output: "))
    output = line.lstrip("// Output: ").strip()

    # Run script and compare
    result = subprocess.run(["./Emplode", test+".emp"], text=True, capture_output=True)
    stdout = result.stdout.strip()
    if stdout == output:
        print("\033[32mPassed:", test, "\033[0m")
        success += 1
    else:
        print("\033[31;1mFailed:", test, "\033[0m")
        print('  Expected: \"', output, '"', sep='')
        print('  Found:    \"', stdout, '"', sep='')
        failure += 1

print()
if failure == 0:
    print("Passed: ", success, ", Failed: ", 0, sep='')
else:
    print("Passed: ", success, ", \033[31;1mFailed: ", failure, "\033[0m", sep='')
