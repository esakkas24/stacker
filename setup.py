from setuptools import setup, find_packages

# Function to read the requirements.txt file
def read_requirements():
    with open('requirements.txt') as req_file:
        return req_file.read().splitlines()

setup(
    name="pistacker",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'stacker=stacker.__init__:run_python_command', 
        ],
    },
    install_requires=read_requirements(), 
)