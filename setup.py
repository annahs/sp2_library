try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

config = {
    'description': 'SP2 common code',
    'author': 'Sarah Hanna',
    'url': 'URL to get it at.',
    'download_url': 'Where to download it.',
    'author_email': 'sarah.hanna@gmail.com',
    'version': '0.1',
    'install_requires': ['nose'],
    'packages': ['sp2_library'],
    'scripts': [],
    'name': 'sp2_library'
}

setup(**config)