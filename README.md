# GladuationResearch

### Using MATLAB and Unity.

#### 0°~360° 偏光版を回転し、観測される色。

![0 5mm_53](https://user-images.githubusercontent.com/57553474/82459059-6ece5300-9af2-11ea-8c03-da1ad8b5e935.png)

#### 偏光版を角度固定し、観測角度を 1°~90° で観測した際に発色。

![0 5mm__0-90](https://user-images.githubusercontent.com/57553474/82459081-74c43400-9af2-11ea-9334-38ac771a166f.png)

## Unity で表示するまでの流れ

#### 縦軸0~90横軸0~360

![img](https://user-images.githubusercontent.com/57553474/83354563-77494800-a394-11ea-80f5-d2759bc67dc2.png)

#### Angle90andRotate360.m で氷の厚さを指定し、画像作成(90x360x3) ※Octave で 15 分かかる

#### img_end.m で ↑ の画像を読み込み、image_end.txt へ出力(1x97200)

#### Unity で読み込み
