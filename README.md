# Prey-and-Predator
Цель. Сформировать практические навыки применения дискретных методов решения задачи Коши на базе
конечно-разностных аппроксимаций

Вариант 6:

![image](https://github.com/user-attachments/assets/83f8730d-5ce5-40c1-8137-3555d0707e93)

![image](https://github.com/user-attachments/assets/190f1f85-71d1-4e1b-824d-f225366ddb25)

Оптимальные параметры:

    if (month >= 3 && month <= 4) {
        // Весна: средний рост зайцев, умеренная активность лис
        alpha = 0.3;
        beta = 0.001;
        gamma = 0.15;
        delta = 0.00008;
    } 
    else if (month >= 5 && month <= 8) {
        // Лето: высокий рост зайцев, высокая активность лис
        alpha = 0.5;
        beta = 0.0012;
        gamma = 0.17;
        delta = 0.0001;
    }
    else if (month >= 9 && month <= 10) {
        // Осень: умеренный спад зайцев, активность лис снижается
        alpha = 0.25;
        beta = 0.0009;
        gamma = 0.16;
        delta = 0.00007;
    }
    else {
        // Зима: сильное снижение численности зайцев, низкая активность лис
        alpha = 0.15;
        beta = 0.0005;
        gamma = 0.2;
        delta = 0.00005;
    }

Результаты моделирования системы Лотки-Вольтерра в течение года:

![image](https://github.com/user-attachments/assets/e5b2340c-fe6d-40b6-b010-c1e2bfb6bc0c)

![image](https://github.com/user-attachments/assets/569b4798-066a-41d7-9f86-03639b02151d)

![image](https://github.com/user-attachments/assets/b658d1c5-5b23-481d-81b4-9dc4d6f292ed)


Выводы: добиться стабильности схемы Адамса на протяжении 10 лет не удалось. Через примерно 4 года система дистабилизиурется настолько, что популяция полностью вымирает. 
Схема Рунге-Кутты остается стабильной и на протяжении 10 лет:

![image](https://github.com/user-attachments/assets/a683e71a-e81d-4856-8267-93d66256171c)
