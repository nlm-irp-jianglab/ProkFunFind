import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='ProkFunFind',
                 version='0.2.2',
                 description='Predict the function based on genome sequences',
                 long_description=long_description,
                 author='Xiaofang Jiang',
                 author_email='xiaofang.jiang@nih.gov',
                 packages=[
                     "ProkFunFind",
                     "ProkFunFind.toolkit",
                     "ProkFunFind.read",
                     "ProkFunFind.report",
                     "ProkFunFind.examine",
                     "ProkFunFind.detect.blast_search",
                     "ProkFunFind.detect.ipr_search",
                     "ProkFunFind.detect.hmmer_search",
                     "ProkFunFind.detect.kofam_search",
                     "ProkFunFind.cluster.DBSCAN"
                 ],
                 # packages=setuptools.find_packages(),
                 license='MIT',
                 python_requires='>=3.7',
                 scripts=["bin/prokfunfind"],
                 install_requires=[
                     'biopython',
                     'scikit-learn',
                     'configparser',
                     'typing',
                     'argparse',
                     'numpy',
                     'six',
                     'bcbio-gff'
                 ],
                 zip_safe=False)
#                     'importlib',
