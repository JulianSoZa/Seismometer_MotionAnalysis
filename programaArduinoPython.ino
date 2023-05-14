const int analogInPin = A0; // entrada analógica al que está conectado el módulo
int sensorValue;
float voltageValue;

void setup() {
  Serial.begin(9600);// initialize serial communications at 9600 bps
}

void loop() {
    sensorValue = analogRead(analogInPin);
    voltageValue = map(sensorValue, 0, 1023, 0, 5000);
    Serial.println(voltageValue);
}
