from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='CisRegSeeker',  # Package name
    version='1.0',  # Version number
    description='CisRegSeeker is a CLI tool for identifying and analyzing cis-regulatory elements in genomic sequences.',  # Short description
    long_description=open('README.md').read(),  # Full description from README.md
    long_description_content_type='text/markdown',  # Content type for long description
    author='[Your Name]',  # Replace with actual author names
    author_email='[Your Email]',  # Replace with actual email
    url='https://github.com/[your-github]/CisRegSeeker',  # Update with your GitHub repo URL
    packages=find_packages(),  # Automatically find all packages in the project
    entry_points={
        'console_scripts': [
            'CisRegSeeker=cisregseeker.main:main',  # Register CLI command
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10',  # Specify Python version requirement
    install_requires=required,  # Dependencies from requirements.txt
)
