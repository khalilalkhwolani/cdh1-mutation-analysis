"""Deep learning models for CDH1 mutation analysis."""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import json
from pathlib import Path

from ..utils.logger import LoggerMixin
from ..utils.config import Config


class LSTMClassifier(LoggerMixin):
    """LSTM-based classifier for CDH1 mutation analysis."""
    
    def __init__(self, config: Config):
        """
        Initialize LSTM classifier.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.model_config = config.get('deep_learning.model')
        self.model = None
        self.history = None
        
    def build_model(self, input_shape: Tuple[int, ...], num_classes: int) -> keras.Model:
        """
        Build LSTM model architecture.
        
        Args:
            input_shape: Input tensor shape
            num_classes: Number of output classes
            
        Returns:
            Compiled Keras model
        """
        self.logger.info(f"Building LSTM model with input shape {input_shape}")
        
        model = keras.Sequential([
            # Input layer
            layers.Input(shape=input_shape),
            
            # LSTM layers
            layers.LSTM(
                self.model_config.get('lstm_units', 128),
                return_sequences=True,
                dropout=self.model_config.get('dropout_rate', 0.2)
            ),
            layers.LSTM(
                self.model_config.get('lstm_units', 128) // 2,
                dropout=self.model_config.get('dropout_rate', 0.2)
            ),
            
            # Dense layers
            layers.Dense(
                self.model_config.get('dense_units', 64),
                activation='relu'
            ),
            layers.Dropout(self.model_config.get('dropout_rate', 0.2)),
            
            # Output layer
            layers.Dense(
                num_classes,
                activation=self.model_config.get('output_activation', 'softmax')
            )
        ])
        
        # Compile model
        model.compile(
            optimizer=keras.optimizers.Adam(
                learning_rate=self.config.get('deep_learning.training.learning_rate', 0.001)
            ),
            loss='categorical_crossentropy' if num_classes > 2 else 'binary_crossentropy',
            metrics=['accuracy', 'precision', 'recall']
        )
        
        self.model = model
        self.logger.info(f"Model built successfully with {model.count_params()} parameters")
        
        return model
        
    def get_model_summary(self) -> str:
        """
        Get model architecture summary.
        
        Returns:
            Model summary string
        """
        if self.model is None:
            return "Model not built yet"
            
        import io
        import sys
        
        # Capture summary output
        old_stdout = sys.stdout
        sys.stdout = buffer = io.StringIO()
        
        try:
            self.model.summary()
            summary = buffer.getvalue()
        finally:
            sys.stdout = old_stdout
            
        return summary
        
    def save_model(self, filepath: Path) -> None:
        """
        Save trained model.
        
        Args:
            filepath: Path to save model
        """
        if self.model is None:
            raise ValueError("No model to save")
            
        self.model.save(filepath)
        self.logger.info(f"Model saved to {filepath}")
        
        # Save model configuration
        config_path = filepath.parent / f"{filepath.stem}_config.json"
        with open(config_path, 'w') as f:
            json.dump(self.model_config, f, indent=2)
            
    def load_model(self, filepath: Path) -> keras.Model:
        """
        Load trained model.
        
        Args:
            filepath: Path to model file
            
        Returns:
            Loaded Keras model
        """
        self.model = keras.models.load_model(filepath)
        self.logger.info(f"Model loaded from {filepath}")
        return self.model
        
    def predict(self, X: np.ndarray) -> np.ndarray:
        """
        Make predictions on input data.
        
        Args:
            X: Input data
            
        Returns:
            Predictions
        """
        if self.model is None:
            raise ValueError("Model not built or loaded")
            
        predictions = self.model.predict(X)
        self.logger.info(f"Generated predictions for {len(X)} samples")
        
        return predictions
        
    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        """
        Get prediction probabilities.
        
        Args:
            X: Input data
            
        Returns:
            Prediction probabilities
        """
        return self.predict(X)
        
    def predict_classes(self, X: np.ndarray) -> np.ndarray:
        """
        Get predicted classes.
        
        Args:
            X: Input data
            
        Returns:
            Predicted class indices
        """
        probabilities = self.predict(X)
        return np.argmax(probabilities, axis=1)
        
    def evaluate_model(self, X_test: np.ndarray, y_test: np.ndarray) -> Dict[str, float]:
        """
        Evaluate model performance.
        
        Args:
            X_test: Test input data
            y_test: Test target data
            
        Returns:
            Evaluation metrics
        """
        if self.model is None:
            raise ValueError("Model not built or loaded")
            
        results = self.model.evaluate(X_test, y_test, verbose=0)
        
        metrics = {}
        for i, metric_name in enumerate(self.model.metrics_names):
            metrics[metric_name] = float(results[i])
            
        self.logger.info(f"Model evaluation completed: {metrics}")
        return metrics
        
    def get_layer_outputs(self, X: np.ndarray, layer_name: str) -> np.ndarray:
        """
        Get outputs from a specific layer.
        
        Args:
            X: Input data
            layer_name: Name of layer to extract outputs from
            
        Returns:
            Layer outputs
        """
        if self.model is None:
            raise ValueError("Model not built or loaded")
            
        # Create intermediate model
        layer = self.model.get_layer(layer_name)
        intermediate_model = keras.Model(
            inputs=self.model.input,
            outputs=layer.output
        )
        
        outputs = intermediate_model.predict(X)
        return outputs
        
    def visualize_architecture(self, output_path: Optional[Path] = None) -> None:
        """
        Visualize model architecture.
        
        Args:
            output_path: Path to save visualization
        """
        if self.model is None:
            raise ValueError("Model not built or loaded")
            
        try:
            keras.utils.plot_model(
                self.model,
                to_file=output_path or "model_architecture.png",
                show_shapes=True,
                show_layer_names=True,
                rankdir='TB'
            )
            self.logger.info(f"Model architecture saved to {output_path}")
        except ImportError:
            self.logger.warning("pydot not available, cannot create model visualization")
            
    def get_model_config(self) -> Dict[str, Any]:
        """
        Get model configuration.
        
        Returns:
            Model configuration dictionary
        """
        if self.model is None:
            return self.model_config
            
        return {
            'architecture': self.model.get_config(),
            'parameters': self.model_config,
            'trainable_params': self.model.count_params(),
            'layers': len(self.model.layers)
        }


class CNNClassifier(LoggerMixin):
    """CNN-based classifier for sequence analysis."""
    
    def __init__(self, config: Config):
        """
        Initialize CNN classifier.
        
        Args:
            config: Configuration object
        """
        self.config = config
        self.model_config = config.get('deep_learning.model')
        self.model = None
        
    def build_model(self, input_shape: Tuple[int, ...], num_classes: int) -> keras.Model:
        """
        Build CNN model architecture.
        
        Args:
            input_shape: Input tensor shape
            num_classes: Number of output classes
            
        Returns:
            Compiled Keras model
        """
        self.logger.info(f"Building CNN model with input shape {input_shape}")
        
        model = keras.Sequential([
            # Input layer
            layers.Input(shape=input_shape),
            
            # Convolutional layers
            layers.Conv1D(64, 3, activation='relu'),
            layers.MaxPooling1D(2),
            layers.Conv1D(128, 3, activation='relu'),
            layers.MaxPooling1D(2),
            layers.Conv1D(256, 3, activation='relu'),
            layers.GlobalMaxPooling1D(),
            
            # Dense layers
            layers.Dense(128, activation='relu'),
            layers.Dropout(0.5),
            layers.Dense(num_classes, activation='softmax')
        ])
        
        model.compile(
            optimizer='adam',
            loss='categorical_crossentropy',
            metrics=['accuracy']
        )
        
        self.model = model
        return model