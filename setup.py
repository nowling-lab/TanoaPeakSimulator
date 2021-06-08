from setuptools import setup

setup(
    name='Tanoa',
    version='0.1',
    description='Artificial peak generator for testing peak calling software',
    url='https://github.com/nowling-lab/TanoaPeakSimulator',
    author='John Peters',
    author_email='John.Geraldo.Peters@gmail.com',
    license='Apache-2.0',
    packages=['tanoa'],
    zip_safe=False,
    python_requires=">=3.6",
    scripts=["bin/tanoa_generate_peaks"]
)
