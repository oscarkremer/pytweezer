from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
     name='pytweezer',  
     version='0.0.1',
     scripts=['pytweezer/mentor.py', 'pymentor/error.py'] ,
     author="Oscar Schmitt Kremer",
     author_email="oscar.s.kremer@hotmail.com",
     description="short description here",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/oscarkremer/pytweezer",
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