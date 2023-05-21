const int analogInPin = A0; // entrada analógica al que está conectado el módulo
int sensorValue;
float voltageValue;
float vok, vik, vo;

void setup() {
  Serial.begin(19200);// initialize serial communications at 9600 bps
  //vok = 0;
  //vik = 0;
}

void loop() {
    sensorValue = analogRead(analogInPin);
    voltageValue = map(sensorValue, 0, 1023, 0, 5000);
    Serial.println(voltageValue);
    /*float S;
    float a;
    a = 1;
    S = 1;
    S = a*sensorValue + (1-a)*S;
    Serial.println(S);*/
    //Transformada z
    //FILTRO DIGITAL 1ER ORDEN
    /*vo=0.4*vok+0.59*vik; //40hz
    //vo=0.07079*vok+0.9292*vik; //90hz
    vok=vo;
    vik=float(sensorValue);/*

    voltageValue = map(vo, 0, 1023, 0, 5000);
    Serial.println(voltageValue);*/
}
