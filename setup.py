from setuptools import setup, find_packages

setup(
    name="stacker",
    version="1.0.0",
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'stacker=stacker.__init__:run_python_command'
        ],
    },
)