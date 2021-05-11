from setuptools import setup

setup(
    name='Tanoa',
    version='0.1',
    description='Artificial peak generator for testing peak calling software',
    url='https://github.com/nowling-lab/TanoaPeakSimulator',
    author='John Peters',
    author_email='John.Geraldo.Peters@gmail.com',
    license='Apache-2.0',
    packages=['Tanoa'],
    zip_safe=False,
    python_requires=">=3.6",
    #install_requires = ["secrets", "random"]. These are default?
    scrips=["bin/tanoa_generate_peaks"]
)
