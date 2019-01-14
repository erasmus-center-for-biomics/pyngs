import setuptools

with open("README.md", "rt") as stream:
    long_description = stream.read()

# prepare the setup
setuptools.setup(
    name="pyngs",
    version="0.2.5",
    author="R.W.W. Brouwer",
    author_email="r.w.w.brouwer@gmail.com",
    description="A package to work with NGS data.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/erasmus-center-for-biomics/pyaln",
    packages=setuptools.find_packages(),
    scripts=[
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Next Generation Sequencing"
    ]
)
