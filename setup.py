import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="eaglepy-jmackereth", # Replace with your own username
    version="0.0.1",
    author="J. Ted Mackereth",
    author_email="j.e.mackereth@bham.ac.uk",
    description="Python based tools for accessing and using the EAGLE simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jmackereth/eaglepy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
