from ms2rescore.feature_generators.base import FeatureGeneratorException


class Ms2pipIncorrectModelException(FeatureGeneratorException):

    def __init__(self, message: str, model: str):
        super().__init__(f"Error: {message}, for model {model}")
