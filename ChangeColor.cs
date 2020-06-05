using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;
using UnityEngine;
using UnityEngine.UI;

public class ChangeColor : MonoBehaviour {

    private GameObject mainCamera;
    private GameObject polarizer;
    private int angle, color_red, color_green, color_blue;
    private string[, , ] rgb = new string[90, 360, 3];

    void Start () {
        mainCamera = Camera.main.gameObject;
        polarizer = GameObject.FindGameObjectWithTag ("Polarizer");

        //fileの読み込み -> 0.5mの画像
        string fileName = "C:/Users/kona/image_end.txt";
        var arrText = new List<string> ();
        StreamReader objReader = new StreamReader (fileName);

        string sLine = ""; // ←一時格納用

        // ファイルから1行ずつ読み込む
        while (sLine != null) {
            sLine = objReader.ReadLine ();
            if (sLine != null)
                arrText.Add (sLine); // addメソッドで追記
        }
        objReader.Close ();

        // --- ArrayListから2次元配列に格納する．---
        int cnt = 0;
        string[, ] array = new string[1, 97200]; // 最初に一列で保存するための配列
        foreach (string sOutput in arrText) {
            string[] temp_line = sOutput.Split (' ');
            foreach (string value in temp_line) {
                if (cnt < 97200) {
                    array[0, cnt] = value;
                    cnt++;
                }
            }
        }

        int column, row, color, temporary;
        for (int i = 0; i < 97200; i++) {
            temporary = i % 270;
            column = temporary / 3;
            row = i / 270;
            color = i % 3;

            rgb[column, row, color] = array[0, i];
        }

        Debug.Log (rgb[53, 0, 0]);
    }

    void Update () {
        angle = (int) mainCamera.transform.localEulerAngles.x;

        if (0 > angle || angle > 90) {
            Color color = new Color32 (255, 255, 255, 0);
            GetComponent<Renderer> ().material.color = color;
        } else {
            //      rgb[観測角度(0~90°), 偏光版角度(0~360°), 色(0~2)]
            color_red = int.Parse (rgb[angle, 0, 0]);
            color_green = int.Parse (rgb[angle, 0, 1]);
            color_blue = int.Parse (rgb[angle, 0, 2]);

            Color color = new Color32 ((byte) color_red, (byte) color_green, (byte) color_blue, 1);
            GetComponent<Renderer> ().material.color = color;
        }

    }
}