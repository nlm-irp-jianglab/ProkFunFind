import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='GutFunFind',
                 version='0.2.0',
                 description='Predict the function based on genome sequences',
                 long_description=long_description,
                 author='Xiaofang Jiang',
                 author_email='xiaofang.jiang@nih.gov',
                 packages=[
                     "GutFunFind",
                     "GutFunFind.toolkit",
                     "GutFunFind.read",
                     "GutFunFind.report",
                     "GutFunFind.examine",
                     "GutFunFind.detect.blast_search",
                     "GutFunFind.detect.ipr_search",
                     "GutFunFind.detect.hmmer_search",
                     "GutFunFind.cluster.DBSCAN"
                     ],
                 #packages=setuptools.find_packages(),
                 license='MIT',
                 python_requires='>=3.7',
                 scripts=["bin/run_GutFunFind.py"],
                 install_requires=[
                     'biopython',
                     'scikit-learn',
                     'configparser',
                     'typing',
                     'importlib',
                     'argparse',
                     'numpy',
                     'six'
                 ],
                 zip_safe=False)
