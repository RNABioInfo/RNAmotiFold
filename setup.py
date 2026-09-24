from setuptools import setup,find_packages

setup(name="RNAmotiFold",
      version="1.0.0",
      description="RNA 3D Motif prediction software",
      url="https://github.com/RNABioInfo/RNAmotiFold",
      download_url="https://github.com/RNABioInfo/RNAmotiFold.git",
      author="Marius Sebeke",
      author_email="marius.sebeke@ibmg.uni-stuttgart.de",
      packages=find_packages(),
      include_package_data=True,
      install_requires=["bio>=1.6.2","requests>=2.32.4"],
      entry_points= {"console_scripts":["RNAmotiFold=RNAmotiFold.py"]},
      zip_safe=False)