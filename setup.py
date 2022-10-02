from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
     name='pymentor',  
     version='0.1.1',
     scripts=['pymentor/mentor.py', 'pymentor/error.py'] ,
     author="Oscar Schmitt Kremer",
     author_email="oscar.s.kremer@hotmail.com",
     description="A package to implement a computational model of Mentor didactic robotic arm.",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/oscarkremer/pymentor",
     packages=find_packages(),
     classifiers=[
         "Programming Language :: Python :: 3",
         "Operating System :: OS Independent",
     ],
    install_requires=[
        'numpy',
        'pytest'
    ]
)