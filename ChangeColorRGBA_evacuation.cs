using System;
using System.IO;
using System.Text;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ChangeColorRGBA : MonoBehaviour {

    private GameObject mainCamera;
    public int angle, color_red, color_green, color_blue, colorrrr;
    public string tmp_r, tmp_g, tmp_b;
    public string[,,] rgb = new string[1, 90, 3];

    private string[] strR;
    // private string[] strG;
    // private string[] strB;

    void Start (){
        mainCamera = Camera.main.gameObject;

        // string fileR = "C:/Users/kona/image_red.txt";
        // string fileG = "C:/Users/kona/image_green.txt";
        // string fileB = "C:/Users/kona/image_blue.txt";

        //  //string fileName = "C:/Users/kona/image.txt";
        //  //StreamReader sr = new StreamReader(fileName);

        // StreamReader srR = new StreamReader(fileR);
        // StreamReader srG = new StreamReader(fileG);
        // StreamReader srB = new StreamReader(fileB);

        //while (srR.Peek() >= 0)
        //{
        //     strR = srR.ReadLine().Split(' ');
        //     strG = srG.ReadLine().Split(' ');
        //     strB = srB.ReadLine().Split(' ');
        //}
        

        // srR.Close();
        // srG.Close();
        // srB.Close();
        string fileName = "C:/Users/kona/image_imp.txt";
            var arrText = new List<string>();
            StreamReader objReader = new StreamReader(fileName);

            string sLine = ""; // ←　一時格納用

            // ファイルから1行ずつ読み込む
            while (sLine != null)
            {
                sLine = objReader.ReadLine();
                if (sLine != null)
                    arrText.Add(sLine); // addメソッドで追記
            }
            objReader.Close();

            // --- ArrayListから2次元配列に格納する．---
            int cnt = 0;
            string[,] array = new string[1, 270]; // 最初に一列で保存するための配列
            foreach (string sOutput in arrText)
            {
                string[] temp_line = sOutput.Split(' ');
                foreach(string value in temp_line)
                {
                    if(cnt < 270)
                    {
                        array[0, cnt] = value;
                        cnt++;
                    }
                }
            }

            int column = 0, row = 0, c = 0;
            
            for(int i = 0; i < 270; i++)
            {
                column = i / 3;
                c = i % 3;
                rgb[row, column, c] = array[0, i];
            }
            string tmp = rgb[0, 89, 2];

            // string str = "";
            // for(int i = 0; i < 90; i++){
            //     str = i.ToString("000");
            //     Debug.Log(str + " : " + rgb[0, i, 1]);
            // }
    }

    void Update () {
        angle = (int)mainCamera.transform.localEulerAngles.x;

        if (0 > angle || angle > 90) {
            Color color = new Color32(255, 255, 255, 0);
            GetComponent<Renderer>().material.color = color;
        } else {
            tmp_r = rgb[0, angle, 0];
            tmp_g = rgb[0, angle, 1];
            tmp_b = rgb[0, angle, 2];

            color_red = int.Parse(tmp_r);
            color_green = int.Parse(tmp_g);
            color_blue = int.Parse(tmp_b);

            //color_red = int.Parse(strR[angle]);
            // color_green = int.Parse(strG[angle]);
            // color_blue = int.Parse(strB[angle]);

            Color color = new Color32((byte)color_red, (byte)color_green, (byte)color_blue, 1);
            GetComponent<Renderer>().material.color = color;
        }
        
    }
}