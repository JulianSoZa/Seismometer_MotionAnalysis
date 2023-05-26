const int analogInPin = A0; // entrada analógica al que está conectado el módulo
int sensorValue;
float voltageValue;

void setup() {
  Serial.begin(19200);
}

void loop() {
    sensorValue = analogRead(analogInPin);
    voltageValue = map(sensorValue, 0, 1023, 0, 5000);
    Serial.println(voltageValue);
}
