import os

from setuptools import setup, find_packages


def readme():
    """
    Utility function to read the README file.
    Used for the long_description.  It's nice, because now 1) we have a top level
    README file and 2) it's easier to type in the README file than to put a raw
    string in below ...
    :return: String
    """
    return open(os.path.join(os.path.dirname(__file__), 'README.md')).read()


setup(
    name='dme_exploration',
    version='0.1.0',
    url='',
    license='',
    author='Jesper Vang/flight505',
    author_email='jesper_vang@me.com',
    description='A streamlit application for identifying patients with DME',
    python_requires='>=3',
    long_description=readme(),
)
