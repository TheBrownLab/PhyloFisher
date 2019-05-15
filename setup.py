import setuptools

setuptools.setup(
    name="msphylo",
    version="0.0.1",
    author="Matt, Alex, Tom, David",
    author_email="matt@example.com",
    description = "blabla",
    packages=['msphylo', 'msphylo_data'],
    package_dir={'msphylo': 'msphylo', 'msphylo_data': 'msphylo/msphylo_data'},
    scripts = ['msphylo/msphylo.py'],
    install_requires=[
          'markdown',
    include_package_data=True

)