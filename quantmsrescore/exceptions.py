from ms2rescore.feature_generators.base import FeatureGeneratorException


class Ms2pipIncorrectModelException(FeatureGeneratorException):

    def __init__(self, message: str, model: str):
        super().__init__(f"Error: {message}, for model {model}")


class MS3NotSupportedException(Exception):

    def __init__(self, message: str):
        super().__init__(f"Error: {message}")


class MzMLNotUnixException(Exception):

    def __init__(self, message: str):
        super().__init__(f"Error: {message}")
